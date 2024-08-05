!  wave speed (eigenvalue) estimation by Goutiere et al. 2008 ASCE

subroutine propagation(sl,sc,sr,u,v,h,g,snm,spec,dd,tsc,poro,hmin,mf)
    implicit none
    double precision, intent(in)  :: u, v, h, g, snm, spec, dd, poro, tsc, hmin, mf
    double precision, intent(out) :: sl, sc, sr

    double precision :: cc, cf, vv, tt, aa

    cc = dsqrt(g*h)
    vv = dsqrt(u**2.d0+v**2.d0)

    if( vv>1e-5 .and. h>hmin ) then
        cf = g*snm**2.d0/h**(1.d0/3.d0)
        tt = cf*vv**2.d0/(spec*g*dd)                    ! Shields number
    else
        cf = 0.d0
        tt = 0.d0
    end if

    if( tt>tsc ) then
        
        !  Meyer, Peter and Muller
        
        !            aa = mf*4.d0/(1.d0-poro)*(tt-tsc)**0.5d0*cc**2.d0*(3.d0*cf*dsqrt(dd/spec/g)*u**3.d0/h/vv   &
        !                 +(tt-tsc)*dsqrt(spec*g*dd**3.d0)*u*v**2.d0/h/(u**2.d0+v**2.d0)**1.5d0)

        !  Ashida and Michiue
!        aa = mf*17.d0/(1.d0-poro)*(dsqrt(tt)-dsqrt(tsc))*cc**2.d0*((3.d0*dsqrt(tt)+dsqrt(tsc))/dsqrt(tt)*cf*dsqrt(dd/spec/g)*u**3.d0/h/vv   &
!                         +(tt-tsc)*dsqrt(spec*g*dd**3.d0)*u*v**2.d0/h/(u**2.d0+v**2.d0)**1.5d0)
        aa = mf*17.d0/(1.d0-poro)*(dsqrt(tt)-dsqrt(tsc))*cc**2.d0*((3.d0*dsqrt(tt)+dsqrt(tsc))/dsqrt(tt)*cf*dsqrt(dd/spec/g)*u*(7.d0*u**2.d0+v**2.d0)/(6.d0*h*vv)   &
                         +(tt-tsc)*dsqrt(spec*g*dd**3.d0)*u*v**2.d0/h/(u**2.d0+v**2.d0)**1.5d0)
                         
        if( u>=0.d0 ) then
            sr = u+cc
            sc = 0.5d0*(u-cc+dsqrt((u-cc)**2.d0+4.d0*aa/(u+cc)))
            sl = 0.5d0*(u-cc-dsqrt((u-cc)**2.d0+4.d0*aa/(u+cc)))
        else
            sr = 0.5d0*(u+cc+dsqrt((u+cc)**2.d0+4.d0*aa/(u-cc)))
            sc = 0.5d0*(u+cc-dsqrt((u+cc)**2.d0+4.d0*aa/(u-cc)))
            sl = u-cc
        end if
    else
        sr = u+cc
        sc = 0.d0
        sl = u-cc
    end if


end subroutine

 ! flux limiter for TVD

subroutine phical(phi,c,dql,dqc,dqr)
    implicit none
    double precision, intent(in)  :: dql,dqc,dqr,c
    double precision, intent(out) :: phi

    double precision :: r, b

    if (c>0) then
        r = dql/dqc
    else
        r = dqr/dqc
    end if
            
!    if (r<0.) then
!        phi = 1.
!    else
!!        phi = 1.-(1.-dabs(c))*r*(1.+r)/(1.+r**2.)  ! vanalbada
!        phi = 1.d0-(1.d0-dabs(c))*2.d0*r/(1.d0+r)   ! vanleer
!    end if

    ! stable  minbee
            
!    if (r<0.) then
!        phi = 1.d0
!    else if ( r>0 .and. r<1. ) then
!        phi = 1.d0-(1.d0-dabs(c))*r
!    else
!        phi = dabs(c)
!    end if
            
     ! superrbee

    if (r<=0.d0) then
        phi = 1.d0
    else if (r>0.d0 .and. r<=0.5d0) then
        phi = 1.d0-2.d0*(1.-dabs(c))*r
    else if (r>0.5d0 .and. r<=1.d0) then
        phi = dabs(c)
    else if (r>1.d0 .and. r<=1.5d0) then
        phi = 1.-(1.-dabs(c))*r
    else
        phi = 2.*dabs(c)-1.d0
    end if

         ! koren

!    if (r<=0.d0) then
!        phi = 1.d0
!    else if (r>0.d0 .and. r<=0.25d0) then
!        phi = 1.d0-2.d0*(1.-dabs(c))*r
!    else if (r>0.25d0 .and. r<=2.5d0) then
!        phi = 1.d0-(1.d0+2.d0*r)*(1.-dabs(c))/3.d0
!    else
!        phi = 2.*dabs(c)-1.d0
!    end if

    ! TCDF

!    if( r<=0.d0 ) then
!        b = r*(1.d0+r)/(1.d0+r*r)
!    else if ( r>0.d0 .and. r<=0.5d0 ) then
!        b = r*r*r-2.d0*r*r+2.d0*r
!    else if ( r>0.5d0 .and. r<=2.d0 ) then
!        b = 0.75d0*r+0.25d0
!    else
!        b = (2.d0*r*r-2.d0*r-2.25d0)/(r*r-r-1.d0)
!    end if

!    phi = 1.d0-(1.d0-dabs(c))*b

end subroutine

subroutine fluxcalx(ffx,ffy,qbx1,qbx2,u,v,h,z,qx,qy,snm,spec,diam,mu_s,tsc,poro,hmin,g,mf,rdx,rdy,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, spec, diam, mu_s, mf, tsc, hmin, g, rdx, rdy, poro
    double precision, dimension(0:nx,0:ny), intent(in) :: h, z, qx, qy
    double precision, dimension(0:nx,0:ny), intent(out) :: ffx, ffy, qbx1, qbx2, u, v

    integer :: i,j
    double precision :: vv, tau, qbs, qbn, dzdx, dzdy, dzdn, uc, vc, hc

          ! calculation of fluxes at grid point

!$omp do private(i,j,vv,tau,qbs)
        do j=0,ny
            do i=0,nx
                if ( h(i,j)>hmin ) then
                    u(i,j) = qx(i,j)/h(i,j)
                    v(i,j) = qy(i,j)/h(i,j)
                    
                    ffx(i,j) = qx(i,j)**2.d0/h(i,j)+0.5d0*g*h(i,j)**2.d0
                    ffy(i,j) = qx(i,j)*qy(i,j)/h(i,j)

                    vv = dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)
                    tau = snm**2.d0*vv**2.d0/(spec*diam*h(i,j)**(1.d0/3.d0))

                    if( tau>tsc .and. vv>1e-5 ) then
!                        qbs = mf*4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                        qbs = mf*17.d0*(tau-tsc)*(dsqrt(tau)-dsqrt(tsc))*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                        qbx1(i,j) = u(i,j)/vv*qbs
                    else
                        qbx1(i,j) = 0.d0
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

            ! calculate bedload flux corresponding to the local slope effect at i+1/2 as source term for Exner eq. 
        
!$omp do private(i,j,uc,vc,hc,vv,tau,dzdx,dzdy,dzdn,qbs,qbn)
        do j=1,ny-1
            do i=0,nx-1
                uc = (u(i,j)+u(i+1,j))*0.5d0
                vc = (v(i,j)+v(i+1,j))*0.5d0
                hc = (h(i,j)+h(i+1,j))*0.5d0

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
!                    qbs = mf*4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                    qbs = mf*17.d0*(tau-tsc)*(dsqrt(tau)-dsqrt(tsc))*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
                    qbx2(i,j) = -vc/vv*qbn
                else
                    qbx2(i,j) = 0.d0
                end if

            end do
        end do

end subroutine

subroutine fluxcaly(ffx,ffy,qby1,qby2,u,v,h,z,qx,qy,snm,spec,diam,mu_s,tsc,poro,hmin,g,mf,rdx,rdy,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, spec, diam, mu_s, mf, tsc, hmin, g, rdx, rdy, poro
    double precision, dimension(0:nx,0:ny), intent(in) :: h, z, qx, qy
    double precision, dimension(0:nx,0:ny), intent(out) :: ffx, ffy, qby1, qby2, u, v

    integer :: i, j
    double precision :: vv, tau, qbs, qbn, dzdx, dzdy, dzdn, uc, vc, hc
        
!$omp do private(i,j,vv,tau,qbs)
        do j=0,ny
            do i=0,nx
                if ( h(i,j)>hmin ) then
                    u(i,j) = qx(i,j)/h(i,j)
                    v(i,j) = qy(i,j)/h(i,j)
                    
                    ffx(i,j) = qx(i,j)*qy(i,j)/h(i,j)
                    ffy(i,j) = qy(i,j)**2.d0/h(i,j)+0.5d0*g*h(i,j)**2.d0
    
                    vv = dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)
                    tau = snm**2.d0*vv**2.d0/(spec*diam*h(i,j)**(1.d0/3.d0))

                    if( tau>tsc .and. vv>1e-5 ) then
!                        qbs = mf*4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                        qbs = mf*17.d0*(tau-tsc)*(dsqrt(tau)-dsqrt(tsc))*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                        qby1(i,j) = v(i,j)/vv*qbs
                    else
                        qby1(i,j) = 0.d0
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
!                    qbs = mf*4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                    qbs = mf*17.d0*(tau-tsc)*(dsqrt(tau)-dsqrt(tsc))*dsqrt(spec*g*diam**3.d0)/(1.d0-poro)
                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
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
        end do

end subroutine

 ! HLL solver for shallow water equation in x-direction

subroutine hllx1storder(f,ff,sl,sr,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sr, u
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j

     ! sl>0 and sr<0 will not take place if bedload transport rate is non-zero

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

subroutine hllxz1storder(f,ff,sl,sc,sr,u,q,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, q
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: cl, cr, ust

!$omp do private(i,j,ust,cl,cr)
    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i+1,j)
            else
                ust = (u(i,j)+u(i+1,j))*0.5d0

                if( dabs(ust)<1e-5 ) then
                    f(i,j) = 0.d0
                else 
                    if (ust>1e-5 ) then
                        cl = sl(i,j) 
                        cr = sc(i,j)
                    else
                        cl = sc(i,j) 
                        cr = sr(i,j)
                    end if

                    f(i,j) = (ff(i,j)*cr-ff(i+1,j)*cl+cr*cl*(q(i+1,j)-q(i,j)))/(cr-cl)
                end if

            end if
        end do
    end do

end subroutine

 ! HLL solver for Exner equation in x-direction

subroutine hllzx1storder(f,ff,sl,sc,sr,z,u,v,snm,dx,qx,h,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, dx
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, v, z, qx, h
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: ust, vst, hst, cl, cr
    
     ! sl>0 and sr<0 will not take place if bedload transport rate is non-zero

!$omp do private(i,j,ust,vst,hst,cl,cr)
    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i+1,j)
            else
                ust = (u(i,j)+u(i+1,j))*0.5d0
                vst = (v(i,j)+v(i+1,j))*0.5d0
                hst = (h(i,j)+h(i+1,j))*0.5d0

                if( dabs(ust)<1e-5 ) then
                    f(i,j) = 0.d0
                else if (ust>1e-5 ) then
                    cl = sl(i,j) 
                    cr = sc(i,j)       

                    f(i,j) = (ff(i,j)*cr-ff(i+1,j)*cl+cr*cl*(z(i+1,j)-z(i,j)))/(cr-cl)
                else
                    cl = sc(i,j) 
                    cr = sr(i,j)

                    f(i,j) = (ff(i,j)*cr-ff(i+1,j)*cl+cr*cl*(z(i+1,j)-z(i,j)))/(cr-cl)
                end if

            end if
        end do
    end do

end subroutine

 ! HLL solver for shallow water equation in y-direction

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

subroutine hllyz1storder(f,ff,sl,sc,sr,u,q,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, q
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: cl, cr, ust

!$omp do private(i,j,ust,cl,cr)
    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i,j+1)
            else
                ust = (u(i,j)+u(i,j+1))*0.5d0

                if( dabs(ust)<1e-5 ) then
                    f(i,j) = 0.d0
                else 
                    if (ust>1e-5 ) then
                        cl = sl(i,j) 
                        cr = sc(i,j)
                    else
                        cl = sc(i,j) 
                        cr = sr(i,j)
                    end if

                    f(i,j) = (ff(i,j)*cr-ff(i,j+1)*cl+cr*cl*(q(i,j+1)-q(i,j)))/(cr-cl)
                end if

            end if
        end do
    end do

end subroutine

 ! HLL solver for Exner equation in x-direction

subroutine hllzy1storder(f,ff,sl,sc,sr,z,u,v,h,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, v, h, z
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: ust, vst, hst, cl, cr

!$omp do private(i,j,ust,vst,hst,cl,cr)
    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i,j+1)
            else
                ust = (u(i,j)+u(i,j+1))*0.5d0
                vst = (v(i,j)+v(i,j+1))*0.5d0
                hst = (h(i,j)+h(i,j+1))*0.5d0

                if( dabs(vst)<1e-5 ) then
                    f(i,j) = 0.d0
                else if (vst>1e-5 ) then
                    cl = sl(i,j) 
                    cr = sc(i,j)
                    
                    f(i,j) = (ff(i,j)*cr-ff(i,j+1)*cl+cr*cl*(z(i,j+1)-z(i,j)))/(cr-cl)
                else
                    cl = sc(i,j) 
                    cr = sr(i,j)

                    f(i,j) = (ff(i,j)*cr-ff(i,j+1)*cl+cr*cl*(z(i,j+1)-z(i,j)))/(cr-cl)
                end if

            end if
        end do
    end do

end subroutine

subroutine waftvdflux(ff,fsta,fl,fr,sl,sr,dql,dqc,dqr,dx,dt)
    implicit none
    double precision, intent(in) :: fsta, fl, fr, sl, sr, dql, dqc, dqr, dx, dt
    double precision, intent(out) :: ff

    double precision :: cl, cr, ull, urr, phi

    if( sl>0.d0 .or. sr<0. ) then
        ff = fsta
    else
        cl = sl*dt/dx
        call phical(phi,cl,dql,dqc,dqr)
        ull = sign(1.d0,cl)*0.5d0*phi*(fsta-fl)
        
        cr = sr*dt/dx
        call phical(phi,cr,dql,dqc,dqr)
        urr = sign(1.d0,cr)*0.5d0*phi*(fr-fsta)

        ff = 0.5d0*(fl+fr)-ull-urr
    end if

end subroutine

subroutine waftvdqb(ff,fsta,fl,fr,sl,sc,sr,dql,dqc,dqr,ul,ur,dx,dt)
    implicit none
    double precision, intent(in) :: fsta, fl, fr, sl, sc, sr, dql, dqc, dqr, ul, ur, dx, dt
    double precision, intent(out) :: ff

    double precision :: cl, cr, ull, urr, phi, ust

    if( sl>0.d0 .or. sr<0. ) then
        ff = fsta
    else
        ust = 0.5d0*(ul+ur)

        if( ust>0.d0 ) then
            cl = sl*dt/dx
            cr = sc*dt/dx
        else
            cl = sc*dt/dx
            cr = sr*dt/dx
        end if

        call phical(phi,cl,dql,dqc,dqr)
        ull = sign(1.d0,cl)*0.5d0*phi*(fsta-fl)
        
        call phical(phi,cr,dql,dqc,dqr)
        urr = sign(1.d0,cr)*0.5d0*phi*(fr-fsta)

        ff = 0.5d0*(fl+fr)-ull-urr
    end if

end subroutine

 ! TVD version of WAF method for shallow water model in x-direction

subroutine tvdwafx(ff,fsta,f,dq,sl,sr,dx,dt,nx,ny)
        implicit none
        integer, intent(in) :: nx, ny
        double precision, intent(in) :: dx, dt
        double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sr,dq
        double precision, dimension(0:nx,0:ny), intent(out) :: ff
    
        integer :: i, j
        double precision :: cr, cl, ull, urr, phi

         !  the case sl >0 or sr<0 uses 1st order HLL flux for boundary
    
!$omp single
        do j=1,ny-1
            i=0             ! for upstream boudary dq(i-1)=dp(i)
            call waftvdflux(ff(i,j),fsta(i,j),f(i,j),f(i+1,j),sl(i,j),sr(i,j),dq(i,j),dq(i,j),dq(i+1,j),dx,dt)

            i=nx-1             ! for downstream boudary dq(i+1)=dp(i)
            call waftvdflux(ff(i,j),fsta(i,j),f(i,j),f(i+1,j),sl(i,j),sr(i,j),dq(i-1,j),dq(i,j),dq(i,j),dx,dt)

        end do
!$omp end single

!$omp do private(i,j,cl,ull,cr,urr,phi)
        do j=1,ny-1
            do i=1,nx-2
                call waftvdflux(ff(i,j),fsta(i,j),f(i,j),f(i+1,j),sl(i,j),sr(i,j),dq(i-1,j),dq(i,j),dq(i+1,j),dx,dt)
            end do
        end do

end subroutine

subroutine tvdwafx_qy(ff,qy,h,f1,dx,dt,nx,ny,hmin,epsilon)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dx, dt, hmin, epsilon
    double precision, dimension(0:nx,0:ny), intent(in) :: qy, f1, h
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j
    double precision :: uc, hc, dqc, dql, dqr, wdq, cc, phi, ucc

!$omp do private(i,j,uc,hc,dqc,dql,dqr,cc,phi,ucc,wdq)
        do j=1,ny-1
            do i=0,nx-1

                hc = (h(i,j)+h(i+1,j))*0.5d0
                if( hc>hmin ) then
                    uc = f1(i,j)/hc

!                    dqc = h(i+1,j)-h(i,j)
                    wdq = qy(i+1,j)-qy(i,j)
                    dqc = wdq
                    if( dabs(wdq)<epsilon ) dqc = epsilon*sign(1.d0,wdq)
                    
                    if( i==0 ) then
                        dql = epsilon
                    else
!                        dql = h(i,j)-h(i-1,j)
                        wdq = qy(i,j)-qy(i-1,j)
                        dql = wdq
                        if( dabs(wdq)<epsilon ) dql = epsilon*sign(1.d0,wdq)
                    end if

                    if( i==nx-1 ) then
                        dqr = epsilon
                    else
!                        dqr = h(i+2,j)-h(i+1,j)
                        wdq = qy(i+2,j)-qy(i+1,j)
                        dqr = wdq
                        if( dabs(wdq)<epsilon ) dqr = epsilon*sign(1.d0,wdq)
                    end if

                    cc = uc*dt/dx
                    call phical(phi,cc,dql,dqc,dqr)
                    ucc = sign(1.d0,cc)*0.5d0*phi*(qy(i+1,j)-qy(i,j))

                    ff(i,j) = uc*( 0.5d0*(qy(i,j)+qy(i+1,j)) - ucc )

                else
                    ff(i,j) = 0.0d0
                end if

            end do
        end do

end subroutine

 ! TVD version of WAF method for Exneq equation in x-direction

subroutine tvdwafxz(ff,fsta,f,z,dq,sl,sc,sr,u,v,h,dx,dt,g,nx,ny)
    implicit none
        integer, intent(in) :: nx, ny
        double precision, intent(in) :: dx, dt, g
        double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sc, sr, dq, u, v, h, z
        double precision, dimension(0:nx,0:ny), intent(out) :: ff
    
        integer :: i, j
        double precision :: cr, cl, ull, urr, phi, ust, vst, hst, fr


!$omp single
        do j=1,ny-1
            i=0
            call waftvdqb(ff(i,j),fsta(i,j),f(i,j),f(i+1,j),sl(i,j),sc(i,j),sr(i,j),dq(i,j),dq(i,j),dq(i+1,j),u(i,j),u(i+1,j),dx,dt)

            i=nx-1
            call waftvdqb(ff(i,j),fsta(i,j),f(i,j),f(i+1,j),sl(i,j),sc(i,j),sr(i,j),dq(i-1,j),dq(i,j),dq(i,j),u(i,j),u(i+1,j),dx,dt)
        end do
!$omp end single

!$omp do private(i,j,ust,vst,hst,fr,cl,cr,phi,ull,urr)
        do j=1,ny-1
            do i=1,nx-2
                call waftvdqb(ff(i,j),fsta(i,j),f(i,j),f(i+1,j),sl(i,j),sc(i,j),sr(i,j),dq(i-1,j),dq(i,j),dq(i+1,j),u(i,j),u(i+1,j),dx,dt)
            end do
        end do

end subroutine

 ! TVD version of WAF method for shallow water model in y-direction

subroutine tvdwafy(ff,fsta,f,dq,sl,sr,dy,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dy, dt
    double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sr, dq
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j
    double precision :: cr, cl, ull, urr, phi

!$omp single
    do i=1,nx-1
        j=0
        call waftvdflux(ff(i,j),fsta(i,j),f(i,j),f(i,j+1),sl(i,j),sr(i,j),dq(i,j),dq(i,j),dq(i,j+1),dy,dt)

        j=ny-1
        call waftvdflux(ff(i,j),fsta(i,j),f(i,j),f(i,j+1),sl(i,j),sr(i,j),dq(i,j-1),dq(i,j),dq(i,j),dy,dt)
    end do
!$omp end single

!$omp do private(i,j,cl,cr,phi,ull,urr)
    do j=1,ny-2
        do i=1,nx-1
            call waftvdflux(ff(i,j),fsta(i,j),f(i,j),f(i,j+1),sl(i,j),sr(i,j),dq(i,j-1),dq(i,j),dq(i,j+1),dy,dt)
        end do
    end do

end subroutine

subroutine tvdwafy_qx(ff,qx,h,f1,dy,dt,nx,ny,hmin,epsilon)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dy, dt, hmin, epsilon
    double precision, dimension(0:nx,0:ny), intent(in) :: qx, f1, h
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j
    double precision :: vc, hc, dqc, dql, dqr, wdq, cc, phi, ucc

!$omp do private(i,j,hc,vc,dqc,dql,dqr,cc,phi,ucc,wdq)
        do j=0,ny-1
            do i=1,nx-1
                                
                hc = (h(i,j)+h(i,j+1))*0.5d0
                if( hc>hmin ) then
                    vc = f1(i,j)/hc

!                    dqc = h(i,j+1)-h(i,j)
                    wdq = qx(i,j+1)-qx(i,j)
                    dqc = wdq
                    if( dabs(wdq)<epsilon ) dqc = epsilon*sign(1.d0,wdq)
                    
                    if( j==0 ) then
                        dql = epsilon
                    else
!                        dql = h(i,j)-h(i,j-1)
                        wdq = qx(i,j)-qx(i,j-1)
                        dql = wdq
                        if( dabs(wdq)<epsilon ) dql = epsilon*sign(1.d0,wdq)
                    end if

                    if( j==ny-1 ) then
                        dqr = epsilon
                    else
!                        dqr = h(i,j+2)-h(i,j+1)
                        wdq = qx(i,j+2)-qx(i,j+1)
                        dqr = wdq
                        if( dabs(wdq)<epsilon ) dqr = epsilon*sign(1.d0,wdq)
                    end if

                    cc = vc*dt/dy
                    call phical(phi,cc,dql,dqc,dqr)
                    ucc = sign(1.d0,cc)*0.5d0*phi*(qx(i,j+1)-qx(i,j))

                    ff(i,j) = vc*( 0.5d0*(qx(i,j)+qx(i,j+1)) - ucc )

                else
                    ff(i,j) = 0.0d0
                end if
            end do
        end do

end subroutine

 ! TVD version of WAF method for Exneq equation in x-direction

subroutine tvdwafyz(ff,fsta,f,z,dq,sl,sc,sr,u,v,h,dy,dt,g,nx,ny)
implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dy, dt, g
    double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sc, sr, dq, u, v, h, z
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j
    double precision :: cr, cl, ull, urr, phi, ust, vst, hst, fr

!$omp single
    do i=1,nx-1
        j=0
        ff(i,j) = 0.d0

        j=ny-1
        ff(i,j) = 0.d0

    end do
!$omp end single

!$omp do private(i,j,ust,vst,hst,fr,cl,cr,phi,ull,urr)
    do j=1,ny-2
        do i=1,nx-1
            call waftvdqb(ff(i,j),fsta(i,j),f(i,j),f(i,j+1),sl(i,j),sc(i,j),sr(i,j),dq(i,j-1),dq(i,j),dq(i,j+1),v(i,j),v(i,j+1),dy,dt)
        end do
    end do

end subroutine

subroutine dqxcal(dh,q,epsilon,nx,ny)
    implicit none
    integer,intent(in) :: nx, ny
    double precision,intent(in) :: epsilon
    double precision, dimension(0:nx,0:ny), intent(in) :: q
    double precision, dimension(0:nx,0:ny), intent(out) :: dh

    integer :: i, j
    double precision :: wdq

!$omp do private(i,j,wdq)
        do j=0,ny
            do i=0,nx-1
                wdq = (q(i+1,j))-(q(i,j))
                dh(i,j) = wdq
                if( dabs(wdq)<epsilon ) dh(i,j) = epsilon*sign(1.d0,wdq)
            end do
        end do

end subroutine

subroutine dqycal(dh,q,epsilon,nx,ny)
    implicit none
    integer,intent(in) :: nx, ny
    double precision,intent(in) :: epsilon
    double precision, dimension(0:nx,0:ny), intent(in) :: q
    double precision, dimension(0:nx,0:ny), intent(out) :: dh

    integer :: i, j
    double precision :: wdq

!$omp do private(i,j,wdq)
    do j=0,ny-1
        do i=0,nx
            wdq = (q(i,j+1))-(q(i,j))
            dh(i,j) = wdq
            if( dabs(wdq)<epsilon ) dh(i,j) = epsilon*sign(1.d0,wdq)
        end do
    end do

end subroutine

subroutine systemx(h,wh,qx,wqx,qy,wqy,z,wz,f1,f2,f3,f4,qbx,qbtsc,u,v,sr,sc,sl,dz,snm,g,poro,hmin,mf,time,bedtime,rdx,rdy,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, g, poro, hmin, mf, time, bedtime, rdx, rdy, dt
    double precision, dimension(0:nx,0:ny), intent(in) :: wh, wqx, wqy, wz, f1, f2, f3, f4, qbx, u, v, sr, sc, sl, qbtsc
    double precision, dimension(0:nx,0:ny), intent(out) :: h, qx, qy, z, dz

    integer :: i, j
    double precision :: roughness, sigl, sigr, advection, wdz, dzc, dzsr, dzsl
    double precision :: cr, cl

!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                h(i,j) = wh(i,j)-(-f1(i-1,j)+f1(i,j))*rdx*dt
                if ( h(i,j)<=hmin ) h(i,j) = hmin
            end do
        end do
    
!$omp do private(i,j,roughness,sigl,sigr,advection,cr,cl,dzc,dzsr,dzsl)
        do j=1,ny-1
            do i=1,nx-1
                roughness = g*snm**2.d0*dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)/wh(i,j)**(4.d0/3.d0)

                sigl = f2(i-1,j)-(sr(i-1,j))/(sr(i-1,j)-sl(i-1,j))*g*(wh(i,j)+wh(i-1,j))*0.5*(wz(i,j)-wz(i-1,j))
                sigr = f2(i,j)-(sl(i,j))/(sr(i,j)-sl(i,j))*g*(wh(i+1,j)+wh(i,j))*0.5*(wz(i+1,j)-wz(i,j))
        
                advection = (-sigl+sigr)*rdx

                qx(i,j) = (wqx(i,j)+(-advection)*dt)/(1.d0+roughness*dt)
            end do
        end do

!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                qy(i,j) = wqy(i,j)-(-f3(i-1,j)+f3(i,j))*rdx*dt
            end do
        end do
        
        if( time>bedtime ) then

!$omp do private(i,j,wdz)
            do j=1,ny-1
                do i=1,nx-1
                    wdz = -((-f4(i-1,j)+f4(i,j))*rdx+(-qbx(i-1,j)+qbx(i,j))*rdx)*dt
                    
                    if( i<0.95*nx ) then
                        z(i,j) = wz(i,j)+wdz
                        dz(i,j) = wdz
                    else
                        z(i,j) = wz(i,j)
                        dz(i,j) = 0.d0
                    end if
                end do
            end do

        else
!$omp do private(i,j)
            do j=1,ny-1
                do i=1,nx-1
                    z(i,j) = wz(i,j)
                    dz(i,j) = 0.d0
                end do
            end do
        end if

end subroutine

subroutine systemy(h,wh,qx,wqx,qy,wqy,z,wz,f1,f2,f3,f4,qby,qbtsc,u,v,sr,sc,sl,dz,snm,g,poro,hmin,mf,time,bedtime,rdx,rdy,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, g, poro, hmin, mf, time, bedtime, rdx, rdy, dt
    double precision, dimension(0:nx,0:ny), intent(in) :: wh, wqx, wqy, wz, f1, f2, f3, f4, qby, u, v, sr, sc, sl, qbtsc
    double precision, dimension(0:nx,0:ny), intent(out) :: h, qx, qy, z, dz

    integer :: i, j
    double precision :: roughness, sigl, sigr, advection, wdz, dzsl, dzsr
    double precision :: cl, cr, dzc
        
!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                h(i,j) = wh(i,j)-(-f1(i,j-1)+f1(i,j))*rdy*dt
                if ( h(i,j)<=hmin ) h(i,j) = hmin
            end do
        end do
           
!$omp do private(i,j,roughness,sigl,sigr,advection,cr,cl,dzc,dzsr,dzsl)
        do j=1,ny-1
            do i=1,nx-1
                roughness = g*snm**2.d0*dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)/wh(i,j)**(4.d0/3.d0)

                sigl = f3(i,j-1)-(sr(i,j-1))/(sr(i,j-1)-sl(i,j-1))*g*(wh(i,j)+wh(i,j-1))*0.5d0*(wz(i,j)-wz(i,j-1))
                sigr = f3(i,j)-(sl(i,j))/(sr(i,j)-sl(i,j))*g*(wh(i,j+1)+wh(i,j))*0.5d0*(wz(i,j+1)-wz(i,j))

                advection = (-sigl+sigr)*rdy

                qy(i,j) = (wqy(i,j)+(-advection)*dt)/(1.d0+roughness*dt)
            end do
        end do

!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                qx(i,j) = wqx(i,j)-(-f2(i,j-1)+f2(i,j))*rdy*dt
            end do
        end do
                
        if( time>bedtime ) then

!$omp do private(i,j,wdz)
            do j=1,ny-1
                do i=1,nx-1
                    wdz = -((-f4(i,j-1)+f4(i,j))*rdy+(-qby(i,j-1)+qby(i,j))*rdy)*dt

                    if( i<0.95*nx ) then
                        z(i,j) = wz(i,j)+wdz
                        dz(i,j) = wdz
                    else
                        z(i,j) = wz(i,j)
                        dz(i,j) = 0.d0
                    end if
                end do
            end do

        else

!$omp do private(i,j,wdz)
            do j=1,ny-1
                do i=1,nx-1
                    z(i,j) = wz(i,j)
                    dz(i,j) = 0.d0
                end do
            end do

        end if

end subroutine

subroutine boundary(h,qx,qy,z,dz,dis,wid,snm,ib,dy,nx,ny)
    implicit none
    integer, intent(in)  :: nx, ny
    double precision, intent(in) :: dis, wid, dy, snm, ib
    double precision, dimension(0:nx,0:ny), intent(inout) :: h, qx, qy, z, dz

    integer :: i, j, nnn
    double precision :: rp, rr
    
!$omp single
    
    call system_clock(count=nnn)
    call random_seed(put=(/nnn/))

    do j=0,ny
        call random_number(rr)
        rp = 1.d0-(rr-0.5d0)*0.05d0

        h( 0,j) = h(   1,j)
        h(nx,j) = h(nx-1,j)

        qx( 0,j) = dis/wid*rp
        qx(nx,j) = qx(nx-1,j)

        qy( 0,j) = qy(   1,j)
        qy(nx,j) = qy(nx-1,j)

        dz( 0,j) = 0.d0
        z( 0,j) = z( 0,j)+dz( 0,j)

        dz(nx,j) = 0.d0 !dz(nx-1,j)
        z(nx,j) = z(nx,j)+dz(nx,j)
        
    end do

    do i=0,nx
        h(i, 0) = h(i,   1)
        h(i,ny) = h(i,ny-1)

        qx(i, 0) = qx(i,   1)
        qx(i,ny) = qx(i,ny-1)

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

  ! Iwagaki's critical shields number calculation

subroutine CriticalShieldsIwagaki(ddd,tsci,spec,snu00,g)
    implicit none
    double precision,intent(in)  :: ddd, spec, snu00, g
    double precision,intent(out) :: tsci
    double precision :: rp_c, rst, usci

    rp_c = dsqrt(spec*g)/snu00
    rst = rp_c*ddd**1.5

    if( rst<=2.14 )	usci = 0.14d0*spec*g*ddd
    if( rst>2.14  )	usci = (0.1235d0*spec*g)**(0.78125)*snu00**(0.4375)*ddd**(0.34375)
    if( rst>54.2  )	usci = 0.034d0*spec*g*ddd
    if( rst>162.7 )	usci = (0.01505d0*spec*g)**(1.136364)*snu00**(-0.272727273)*ddd**(1.40909091)
    if( rst>671   )	usci = 0.05d0*spec*g*ddd

    usci = dsqrt(usci)
    tsci = usci**2.d0/(spec*g*ddd)

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
    double precision :: dql, dqc, dqr, ucc, cc, phi, epsilon, wdq
    double precision :: mf, ctm, dtm, sm, dum
    double precision, dimension(:,:), allocatable :: x, y, z, wz, z0, dz, h, wh, u, v, qx, wqx, qy, wqy, sr, sl, sc
    double precision, dimension(:,:), allocatable :: f1, f2, f3, ffx, ffy, qbsta, ffxz, ffyz
    double precision, dimension(:,:), allocatable :: qbx, qby, qbx1, qby1, qbx2, qby2, f4
    double precision, dimension(:,:), allocatable :: qsta, fsta, dh, qbtsc
    character(20) :: fpout, fptxt

    g  = 9.81d0
    nu = 1e-6
    pi = 3.14159d0
    epsilon = 1e-9

    dis     = 0.01035d0     ! water discharge (m3/s)
!    dis = 0.0064
    chlen   = 120.d0        ! channel length (m)
    wid     = 0.9d0         ! channel width (m)
    ib      = 0.0125d0      ! bed slope
!    ib      = 0.005d0      ! bed slope
    spec    = 1.65d0        ! specific weight of sediment in fluid
    diam    = 0.00076d0     ! sediment diameter (m)
    poro    = 0.4d0         ! porosity of bed
    mu_s    = 0.7d0         ! dynamic friction coefficinet of sediment
    hmin    = 0.0001d0      ! minimun water depth (m)
    mf = 1.d0               ! morphological acceleration factor

    call CriticalShieldsIwagaki(diam,tsc,spec,nu,g)
    snm = (2.5d0*diam)**(1.d0/6.d0)/(7.66d0*dsqrt(g))

    nx = 1500               ! number of grid point in downstream direction
    ny = 18                 ! number of grid point in transverse direction
    dx = chlen/dble(nx)
!    dy = wid/dble(ny)
    dy = wid/dble(ny-1)
    dt = 0.02d0
    ct = 0.5d0             ! Courant number
    ctm = ct*0.75d0
    
    rdx = 1.d0/dx
    rdy = 1.d0/dy
    
    tuk   = 60.    !/mf             ! output time interval (sec)
    bedtime = 60.           ! start time for morphological change of bed
    etime = tuk*500.        ! end time of calculation

    iomp = 4                ! number of cores for OpenMP parallelization

        ! allocation of arrays

    allocate( x(0:nx,0:ny), y(0:nx,0:ny), z(0:nx,0:ny), wz(0:nx,0:ny), z0(0:nx,0:ny), dz(0:nx,0:ny), h(0:nx,0:ny), wh(0:nx,0:ny), qbsta(0:nx,0:ny) )
    allocate( u(0:nx,0:ny), v(0:nx,0:ny), qx(0:nx,0:ny), wqx(0:nx,0:ny), qy(0:nx,0:ny), wqy(0:nx,0:ny) )
    allocate( sr(0:nx,0:ny), sl(0:nx,0:ny), sc(0:nx,0:ny), ffx(0:nx,0:ny), ffy(0:nx,0:ny), f1(0:nx,0:ny), f2(0:nx,0:ny), f3(0:nx,0:ny) )
    allocate( qbx(0:nx,0:ny), qby(0:nx,0:ny), qbx1(0:nx,0:ny), qby1(0:nx,0:ny), qbx2(0:nx,0:ny), qby2(0:nx,0:ny), f4(0:nx,0:ny) )
    allocate( qsta(0:nx,0:ny), fsta(0:nx,0:ny), dh(0:nx,0:ny) )
    allocate( ffxz(0:nx,0:ny), ffyz(0:nx,0:ny), qbtsc(0:nx,0:ny) )

        ! grid generation. just straight, rectangler grid

    do j=0,ny
        x(0,j) = 0.d0
        do i=1,nx
            x(i,j) = x(i-1,j)+dx
        end do
    end do
    
    do i=0,nx
        y(i,0) = -dy*0.5d0
        do j=1,ny
            y(i,j) = y(i,j-1)+dy
        end do
    end do
       
    z00  = chlen*ib     ! bed elevation at upstream end

    do j=0,ny
        z(0,j) = z00
        do i=1,nx
            z(i,j) = z(i-1,j)-ib*dx
            h(i,j) = (snm*dis/(wid*ib**0.5))**0.6
            qx(i,j) = dis/wid
            qy(i,j) = 0.d0
        end do
    end do
    
    write(*,*) "Fr= ", h(0.5*nx,0.5*ny)**(2./3.)*ib**0.5/snm/dsqrt(g*h(0.5*nx,0.5*ny))

    do j=0,ny
        do i=0,nx
            z0(i,j) = z(i,j)
        end do
    end do

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
    dz = 0.d0
    
    tt = 0.d0           ! time
    optime = 0.d0
    
	m = 0
    fpout = 'vc000.vtk'
    fptxt = 'tt000.txt'

    call system_clock(cal_t1)	! start time of calculation

!$	call omp_set_num_threads(iomp)

!$omp parallel
    
    do          ! start time loop
    
!$omp do private(i,j)
        do j=0,ny
            do i=0,nx
                h(i,j) = wh(i,j)
                z(i,j) = wz(i,j)
                qx(i,j) = wqx(i,j)
                qy(i,j) = wqy(i,j)
            end do
        end do

            !  x direction 

        call fluxcalx(ffx,ffy,qbx1,qbx2,u,v,wh,wz,wqx,wqy,snm,spec,diam,mu_s,tsc,poro,hmin,g,mf,rdx,rdy,nx,ny)

            ! wave speed calculation at i+1/2

!$omp do private(i,j,wsl1,wsc1,wsr1,wsl2,wsc2,wsr2,uc,vc,hc)
        do j=1,ny-1
            do i=0,nx-1
                call propagation(wsl1, wsc1, wsr1,u(i,j),v(i,j),wh(i,j),g,snm,spec,diam,tsc,poro,hmin,mf)
                call propagation(wsl2, wsc2, wsr2,u(i+1,j),v(i+1,j),wh(i+1,j),g,snm,spec,diam,tsc,poro,hmin,mf)

                sl(i,j) = min(wsl1,wsl2)
                sr(i,j) = max(wsr1,wsr2)
                
                if( dabs(wsc1)>dabs(wsc2) ) then
                    sc(i,j) = wsc1
                else
                    sc(i,j) = wsc2
                end if

            end do
        end do

!$omp single
        dt1 = 999.d0
        do j=1,ny-1
            do i=0,nx-1
                dt1 = min(dabs(dx/sl(i,j)), dabs(dx/sr(i,j)),dt1)
            end do
        end do
!$omp end single

            !  HLL solver for estimating flux at star region

        call hllx1storder(f1,wqx,sl,sr,wh,nx,ny)
        call hllx1storder(f2,ffx,sl,sr,wqx,nx,ny)

!$omp do private(i,j,hc,ust,dqc,dql,dqr,cc,phi,ucc,wdq)
        do j=1,ny-1
            do i=0,nx-1
                ust = 0.5d0*(u(i,j)+u(i+1,j))
                if(ust>0) then
                    f3(i,j) = f1(i,j)*v(i,j)
                else
                    f3(i,j) = f1(i,j)*v(i+1,j)
                end if
            end do
        end do

        call hllzx1storder(f4,qbx1,sl,sc,sr,wz,u,v,snm,dx,wqx,wh,nx,ny)

        call systemx(h,wh,qx,wqx,qy,wqy,z,wz,f1,f2,f3,f4,qbx2,qbtsc,u,v,sr,sc,sl,dz,snm,g,poro,hmin,mf,tt,bedtime,rdx,rdy,dt,nx,ny)

        call boundary(h,qx,qy,z,dz,dis,wid,snm,ib,dy,nx,ny)

            ! y direction for dt

        call fluxcaly(ffx,ffy,qby1,qby2,u,v,h,z,qx,qy,snm,spec,diam,mu_s,tsc,poro,hmin,g,mf,rdx,rdy,nx,ny)
        
!$omp do private(i,j,wsl1,wsc1,wsr1,wsl2,wsc2,wsr2,uc,vc,hc)
        do j=0,ny-1
            do i=1,nx-1
                call propagation(wsl1, wsc1, wsr1,v(i,j),u(i,j),h(i,j),g,snm,spec,diam,tsc,poro,hmin,mf)
                call propagation(wsl2, wsc2, wsr2,v(i,j+1),u(i,j+1),h(i,j+1),g,snm,spec,diam,tsc,poro,hmin,mf)

                sl(i,j) = min(wsl1,wsl2)
                sr(i,j) = max(wsr1,wsr2)
                
                if( dabs(wsc1)>dabs(wsc2) ) then
                    sc(i,j) = wsc1
                else
                    sc(i,j) = wsc2
                end if
            end do
        end do

!$omp single
        dt2 = 999.d0
        do j=0,ny-1
            do i=1,nx-1
                dt2 = min(dabs(dy/sl(i,j)), dabs(dy/sr(i,j)),dt2)
            end do
        end do
!$omp end single

        call hlly1storder(f1,qy,sl,sr,h,nx,ny)
        call hlly1storder(f3,ffy,sl,sr,qy,nx,ny)

!$omp single
        do i=0,nx
            f1(i,0) = 0.d0
            f1(i,ny-1) = 0.d0
        end do
!$omp end single

        
!$omp do private(i,j,hc,ust,dqc,dql,dqr,cc,phi,ucc,wdq)
        do j=0,ny-1
            do i=1,nx-1
                ust = 0.5d0*(v(i,j)+v(i,j+1))
                if(ust>0) then
                    f2(i,j) = f1(i,j)*u(i,j)
                else
                    f2(i,j) = f1(i,j)*u(i,j+1)
                end if
            end do
        end do

        call hllzy1storder(f4,qby1,sl,sc,sr,z,u,v,h,nx,ny)
        
        call systemy(wh,h,wqx,qx,wqy,qy,wz,z,f1,f2,f3,f4,qby2,qbtsc,u,v,sr,sc,sl,dz,snm,g,poro,hmin,mf,tt,bedtime,rdx,rdy,dt,nx,ny)

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

            write(*,'(f12.5,f10.6, f10.3)') tt, dt, mf

            write(fpout(3:5),'(i3)') m
            
            do ii=3,4
                if(fpout(ii:ii) == ' ') fpout(ii:ii) = '0'
            end do

            write(fptxt(3:5),'(i3)') m
            
            do ii=3,4
                if(fptxt(ii:ii) == ' ') fptxt(ii:ii) = '0'
            end do

!            call out2paraview(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
            call out2paraview(x,y,z,h,qbx1,qby1,z0,nx,ny,fpout,tt)
!            call out2txt(x,y,z,h,u,v,z0,nx,ny,fptxt,tt)
            call out2txt(x,y,z,h,qbx1,qby1,z0,nx,ny,fptxt,tt)

            if( m>0 ) optime = optime-tuk
            m = m+1
        
        end if
        
        optime = optime+dt
        tt = tt+dt
!$omp end single

        if( tt>etime ) exit
        
    end do

!$omp end parallel

    call system_clock(cal_t2, t_rate, t_max)
	if ( cal_t2 < cal_t1 ) then
		t_diff = (t_max - cal_t1) + cal_t2 + 1
	else
		t_diff = cal_t2 - cal_t1
	endif
	write(*,*) "Calcuration time",real(t_diff/t_rate),"sec."
	write(*,*) "Calcuration time",real(t_diff/t_rate)/60.,"min."
	write(*,*) "Calcuration time",real(t_diff/t_rate)/3600.,"hour."

end program SWE_HLL