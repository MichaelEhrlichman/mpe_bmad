module linear_aperture_mod

use bmad

implicit none

contains

subroutine linear_aperture(ring,da_config)
  !works OK with dispersion at injection point.  TESTED!
  !assumes all apertures are elliptical
  use custom_dynamic_aperture_mod, only: custom_aperture_scan_struct
  use transfer_map_mod

  implicit none

  type(lat_struct) ring
  type(custom_aperture_scan_struct) da_config

  integer i, j, p
  real(rp) tmat(6,6), t1(6,6)
  real(rp) tvec(6), v1(6)
  real(rp) delta_angle
  real(rp) x1,y1,x2,y2
  real(rp) dE, theta
  real(rp) a,b
  real(rp) x1ds,y1ds,x2ds,y2ds
  real(rp) m,k,pl2,det,Ex1,Ey1,Ex2,Ey2
  real(rp) da2
  real(rp) xa, ya, xb, yb, xpick, ypick
  real(rp) xaoff, xboff, yaoff, yboff
  real(rp) vec_length
  real(rp) offx, offy
  logical mask_x(1:ring%n_ele_track)
  logical mask_y(1:ring%n_ele_track)
  integer, allocatable :: quadrants(:)

  allocate(quadrants(da_config%n_angle))

  delta_angle = (da_config%max_angle - da_config%min_angle)/(da_config%n_angle -1)
  if(delta_angle .ne. delta_angle) delta_angle = 0.0d0  !case of n_angle==1.  In fortran NaN does not equal itself.
  a = ring%ele(1)%value(x1_limit$)
  b = ring%ele(1)%value(y1_limit$)

  mask_x = abs(ring%ele(1:ring%n_ele_track)%value(x1_limit$)) > 0.0001
  mask_y = abs(ring%ele(1:ring%n_ele_track)%value(y1_limit$)) > 0.0001
  da_config%Sx = minval( ring%ele(1:ring%n_ele_track)%value(x1_limit$) / sqrt(ring%ele(1:ring%n_ele_track)%a%beta), mask_x ) * sqrt(ring%ele(1)%a%beta)
  da_config%Sy = minval( ring%ele(1:ring%n_ele_track)%value(y1_limit$) / sqrt(ring%ele(1:ring%n_ele_track)%b%beta), mask_y ) * sqrt(ring%ele(1)%b%beta)

  !initialize to physical aperture at track point
  do i=1,da_config%n_angle
    theta = (i-1)*delta_angle + da_config%min_angle
    if(theta .gt. 0 .and. theta .le. pi/2.0d0) then
      quadrants(i) = 1
    elseif(theta .gt. pi/2.0d0 .and. theta .le. pi) then
      quadrants(i) = 2
    elseif(theta .gt. pi .and. theta .le. 3.0*pi/2.0) then
      quadrants(i) = 3
    elseif(theta .gt. 3.0*pi/2.0 .and. theta .le. 2.0*pi) then
      quadrants(i) = 4
    else
      write(*,*) "theta error. theta = ", theta
    endif

    vec_length = sqrt(a**2 * cos(theta)**2 + a**2 * sin(theta)**2)
    da_config%aperture(i)%x = vec_length * local_cosphi(theta,da_config%Sx,da_config%Sy)
    da_config%aperture(i)%y = vec_length * local_sinphi(theta,da_config%Sx,da_config%Sy)
    da_config%aperture(i)%i_turn = 1  !this routine does not calculate the turn number at which the particle was lost.
  enddo

  do i=1,ring%n_ele_track
    a = ring%ele(i)%value(x1_limit$)  !horizontal axis of ellipse
    b = ring%ele(i)%value(y1_limit$)  !vertical axis of ellipse

    call transfer_matrix_calc (ring, t1, v1, ix1=i, one_turn=.true.)
    call transfer_matrix_calc (ring, tmat, tvec, ix1=0, ix2=i)
    do j=1,100  !simulate 100 turns
      if (j>1) then
        call concat_transfer_mat (t1, v1, tmat, tvec, tmat, tvec)
      endif
      do p=1,da_config%n_angle
        theta = (p-1)*delta_angle + da_config%min_angle  !unscaled angle at injection point

        x1 = local_cosphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(1)  ! cos (scaled theta at injection point)
        y1 = local_sinphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(3)  ! sin (scaled theta at injection point)
        x2 = 0.5*local_cosphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(1)  ! cos (scaled theta at injection point)
        y2 = 0.5*local_sinphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(3)  ! sin (scaled theta at injection point)

        offx = tmat(1,6)*da_config%param%closed_orbit%vec(6) + tvec(1)
        offy = tmat(3,6)*da_config%param%closed_orbit%vec(6) + tvec(3)

        x1ds = tmat(1,1)*x1 + tmat(1,3)*y1 + offx  ! horiz component of (xth,yth) transported downstream
        y1ds = tmat(3,1)*x1 + tmat(3,3)*y1 + offy  ! vertical component of (xth,yth) transported downstream
        x2ds = tmat(1,1)*x2 + tmat(1,3)*y2 + offx  ! horiz component of (xth,yth) transported downstream
        y2ds = tmat(3,1)*x2 + tmat(3,3)*y2 + offy  ! vertical component of (xth,yth) transported downstream

        ! y=mx+k form of downstream vector
        m = (y2ds-y1ds)/(x2ds-x1ds)
        k = y1ds - m*x1ds

        Ex1 = a * (-a*k*m-b*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)
        Ex2 = a * (-a*k*m+b*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)
        Ey1 = b * (b*k-a*m*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)
        Ey2 = b * (b*k+a*m*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)

        det = tmat(1,1)*tmat(3,3) - tmat(1,3)*tmat(3,1)
        xa = (tmat(3,3)*(Ex1-offx) - tmat(1,3)*(Ey1-offy))/det  !aperture projected back to injection point
        ya = (-tmat(3,1)*(Ex1-offx) + tmat(1,1)*(Ey1-offy))/det  !aperture projected back to injection point
        xb = (tmat(3,3)*(Ex2-offx) - tmat(1,3)*(Ey2-offy))/det  !aperture projected back to injection point
        yb = (-tmat(3,1)*(Ex2-offx) + tmat(1,1)*(Ey2-offy))/det  !aperture projected back to injection point

        xaoff = xa - da_config%param%closed_orbit%vec(1)
        xboff = xb - da_config%param%closed_orbit%vec(1)
        yaoff = ya - da_config%param%closed_orbit%vec(3)
        yboff = yb - da_config%param%closed_orbit%vec(3)
        if(quadrants(p) == 1) then
          if( yaoff > 0 .and. xaoff > 0) then
            xpick = xa
            ypick = ya
          elseif( yboff > 0 .and. xboff > 0) then
            xpick = xb
            ypick = yb
          else
            write(*,'(a,i6,4f11.4)') "theta error: ", quadrants(p), xaoff, yaoff, xboff, yboff
            stop
          endif
        elseif(quadrants(p) == 2) then
          if( yaoff > 0 .and. xaoff < 0) then
            xpick = xa
            ypick = ya
          elseif( yboff > 0 .and. xboff < 0) then
            xpick = xb
            ypick = yb
          else
            write(*,'(a,i6,4f11.4)') "theta error: ", quadrants(p), xaoff, yaoff, xboff, yboff
            stop
          endif
        elseif(quadrants(p) == 3) then
          if( yaoff < 0 .and. xaoff < 0) then
            xpick = xa
            ypick = ya
          elseif( yboff < 0 .and. xboff < 0) then
            xpick = xb
            ypick = yb
          else
            write(*,'(a,i6,4f11.4)') "theta error: ", quadrants(p), xaoff, yaoff, xboff, yboff
            stop
          endif
        elseif(quadrants(p) == 4) then
          if( yaoff < 0 .and. xaoff > 0) then
            xpick = xa
            ypick = ya
          elseif( yboff < 0 .and. xboff > 0) then
            xpick = xb
            ypick = yb
          else
            write(*,'(a,i6,4f11.4)') "theta error: ", quadrants(p), xaoff, yaoff, xboff, yboff
            stop
          endif
        endif

        pl2 = (xpick-da_config%param%closed_orbit%vec(1))**2 + (ypick-da_config%param%closed_orbit%vec(3))**2
        da2 = (da_config%aperture(p)%x-da_config%param%closed_orbit%vec(1))**2 + (da_config%aperture(p)%y-da_config%param%closed_orbit%vec(3))**2

        if( pl2 .lt. da2 ) then
          da_config%aperture(p)%x = xpick
          da_config%aperture(p)%y = ypick
        endif
      enddo
    enddo
  enddo

  deallocate(quadrants)
end subroutine

function local_cosphi(th,Sx,Sy) result(x)
  real(rp) th, Sx, Sy, x
  x = Sx * cos(th) / sqrt(Sx**2 * cos(th)**2 + Sy**2 * sin(th)**2 )
end function

function local_sinphi(th,Sx,Sy) result(x)
  real(rp) th, Sx, Sy, x
  x = Sy * sin(th) / sqrt(Sx**2 * cos(th)**2 + Sy**2 * sin(th)**2 )
end function


end module





